import os
import subprocess
from collections import OrderedDict
from collections.abc import Iterator

import wtl.options as wopt
from wtl import cli


def iter_values() -> Iterator[OrderedDict[str, str]]:
    crossing_axes: OrderedDict[str, list[str]] = OrderedDict()
    crossing_axes["xi"] = ["1e-4", "1e-3"]
    for d in wopt.product(crossing_axes):
        yield OrderedDict(**d)


def iter_args(
    rest: list[str], concurrency: int, repeat: int, skip: int
) -> Iterator[list[str]]:
    const = ["tek2", f"-j{concurrency}", *rest]
    const.extend(["-i10000", "-g10000", "-s100000"])
    now = wopt.now()
    for i, v in enumerate(wopt.cycle(iter_values(), repeat)):
        if i < skip:
            continue
        args = wopt.make_args(v)
        label = "_".join([wopt.join(args), now, format(i, "03d")])
        yield const + args + ["--outdir=" + label]


def main() -> None:
    parser = wopt.ArgumentParser()
    parser.add_argument("-p", "--parallel", type=int, default=1)
    (args, rest) = parser.parse_known_args()
    print(f"{os.cpu_count()=}")
    print(f"{args.jobs} jobs * {args.parallel} threads/job")

    it = iter_args(rest, args.parallel, args.repeat, args.skip)
    if cli.dry_run:
        for cmd in it:
            print(" ".join(cmd))
    else:
        fs = [cli.thread_submit(subprocess.run, x, check=True) for x in it]
        cli.wait_raise(fs)
    print("End of " + __file__)


if __name__ == "__main__":
    main()
