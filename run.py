#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
import wtl.options as wopt

program = 'tek'


def iter_values(rest):
    parallel_axes = wopt.OrderedDict()
    parallel_axes['alpha'] = ['0.70', '0.75', '0.80', '0.85']
    parallel_axes['beta'] = [6, 12, 24, 48]
    crossing_axes = wopt.OrderedDict()
    crossing_axes['xi'] = ['1e-5', '1e-4', '1e-3']
    crossing_axes['lambda'] = ['1e-4', '1e-3']
    crossing_axes['nu'] = ['0', '1e-6', '1e-4']
    for di in wopt.parallel(parallel_axes):
        for dj in wopt.product(crossing_axes):
            yield wopt.OrderedDict(**di, **dj)


def iter_args(rest, concurrency, repeat, skip):
    const = [program, '-j{}'.format(concurrency)] + rest
    suffix = '_{}'.format(wopt.now())
    for i, v in enumerate(wopt.cycle(iter_values(rest), repeat)):
        if i < skip:
            continue
        args = wopt.make_args(v)
        label = wopt.join(args) + suffix + '_{:02}'.format(i)
        yield const + args + ['--outdir=' + label]


def main():
    parser = wopt.ArgumentParser()
    parser.add_argument('--skip', type=int, default=0)
    parser.add_argument('-o', '--outdir', default='.stdout')
    parser.add_argument('-r', '--repeat', type=int, default=1)
    (args, rest) = parser.parse_known_args()
    print("cpu_count(): {}".format(wopt.cpu_count()))
    print('{} jobs * {} threads/job'.format(args.jobs, args.parallel))

    wopt.map_async(iter_args(rest, args.parallel, args.repeat, args.skip),
                   args.jobs, args.dry_run, outdir=args.outdir)
    print('End of ' + __file__)


if __name__ == '__main__':
    main()
