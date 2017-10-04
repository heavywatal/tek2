#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
import wtl.options as wopt

program = 'tek'


def iter_values(rest):
    crossing_axes = wopt.OrderedDict()
    crossing_axes['xi'] = ['1e-4', '5e-4']
    crossing_axes['spec'] = ['1e-6', '1e-5']
    for d in wopt.product(crossing_axes):
        yield wopt.OrderedDict(**d)


def iter_args(rest, concurrency, repeat, skip):
    const = [program, '-j{}'.format(concurrency)] + rest
    suffix = '_{}_{}'.format(wopt.now(), wopt.getpid())
    for i, v in enumerate(wopt.cycle(iter_values(rest), repeat)):
        if i < skip:
            continue
        args = wopt.make_args(v)
        label = wopt.join(args) + suffix + '_{:02}'.format(i)
        yield const + args + ['--outdir=' + label]


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--dry-run', action='store_true')
    parser.add_argument('-j', '--jobs', type=int, default=wopt.cpu_count())
    parser.add_argument('--skip', type=int, default=0)
    parser.add_argument('-o', '--outdir', default='.stdout')
    parser.add_argument('-r', '--repeat', type=int, default=1)
    (args, rest) = parser.parse_known_args()

    per_job = wopt.cpu_count() // args.jobs
    if args.jobs > 1:
        per_job += 1
    wopt.map_async(iter_args(rest, per_job, args.repeat, args.skip),
                   args.jobs, args.dry_run, outdir=args.outdir)
    print('End of ' + __file__)


if __name__ == '__main__':
    main()
