import os
import subprocess
from collections import OrderedDict
from collections.abc import Callable, Iterable, Iterator, Mapping
from typing import Any

import wtl.options as wopt
from wtl import cli


def te2fig1() -> Iterator[OrderedDict[str, Any]]:
    crossing_axes: OrderedDict[str, list[str]] = OrderedDict()
    crossing_axes["n"] = ["1000"]
    crossing_axes["xi"] = ["10e-4", "5e-4", "1e-4"]
    crossing_axes["g"] = ["50000"]
    crossing_axes["i"] = ["100"]
    for d in wopt.product(crossing_axes):
        yield OrderedDict(**d)


def te2fig2() -> Iterator[OrderedDict[str, Any]]:
    crossing_axes: OrderedDict[str, list[str]] = OrderedDict()
    crossing_axes["n"] = ["1000"]
    crossing_axes["xi"] = ["10e-4"]
    crossing_axes["g"] = ["50000"]
    crossing_axes["i"] = ["100"]
    for d in wopt.product(crossing_axes):
        yield OrderedDict(**d)


def te2fig4() -> list[dict[str, str]]:
    return [
        {
            "n": "1000",
            "xi": "10e-4",
            "g": "6000",
            "H": "4000",
            "i": "20",
        }
    ]


def te2fig5() -> Iterator[OrderedDict[str, Any]]:
    crossing_axes: OrderedDict[str, list[str]] = OrderedDict()
    crossing_axes["n"] = ["1000"]
    crossing_axes["xi"] = ["10e-4"]
    crossing_axes["coexist"] = ["2"]
    crossing_axes["lower"] = ["6", "9"]
    crossing_axes["upper"] = ["18", "24", "30"]
    crossing_axes["g"] = ["50000"]
    crossing_axes["i"] = ["100"]
    for d in wopt.product(crossing_axes):
        yield OrderedDict(**d)
        yield OrderedDict(**d)


def te2fig6() -> Iterator[OrderedDict[str, Any]]:
    crossing_axes: OrderedDict[str, list[str]] = OrderedDict()
    crossing_axes["r"] = ["1"]
    crossing_axes["n"] = ["1000"]
    crossing_axes["xi"] = ["10e-4"]
    crossing_axes["coexist"] = ["2", "5", "8"]
    crossing_axes["lower"] = ["6", "9"]
    crossing_axes["upper"] = ["18", "24", "30"]
    crossing_axes["g"] = ["50000"]
    crossing_axes["i"] = ["100"]
    for d in wopt.product(crossing_axes):
        yield OrderedDict(**d)


def te1fig2s() -> Iterator[OrderedDict[str, Any]]:
    parallel_axes: OrderedDict[str, list[str | int]] = OrderedDict()
    parallel_axes["alpha"] = ["0.70", "0.75", "0.80", "0.85"]
    parallel_axes["beta"] = [6, 12, 24, 48]
    crossing_axes: OrderedDict[str, list[str]] = OrderedDict()
    crossing_axes["xi"] = ["1e-5", "1e-4", "1e-3"]
    crossing_axes["lambda"] = ["1e-4", "1e-3"]
    crossing_axes["nu"] = ["0", "1e-6", "1e-4"]
    for di in wopt.parallel(parallel_axes):
        for dj in wopt.product(crossing_axes):
            yield OrderedDict(**di, **dj)


def iter_args(
    arg_maker: Callable[..., Iterable[Mapping[str, Any]]],
    rest: list[str],
    concurrency: int,
    repeat: int,
    skip: int,
) -> Iterator[list[str]]:
    const = ["tek2", f"-j{concurrency}", *rest]
    now = wopt.now()
    for i, v in enumerate(wopt.cycle(arg_maker(), repeat)):
        if i < skip:
            continue
        args = wopt.make_args(v)
        label = "_".join([wopt.join(args), now, format(i, "03d")])
        yield const + args + ["--outdir=" + label]


def main() -> None:
    arg_makers = {
        k: v
        for k, v in globals().items()
        if callable(v) and k not in ("main", "iter_args")
    }
    parser = wopt.ArgumentParser()
    parser.add_argument("-p", "--parallel", type=int, default=1)
    parser.add_argument("function", choices=arg_makers)
    (args, rest) = parser.parse_known_args()
    print(f"{os.cpu_count()=}")
    print(f"{args.jobs} jobs * {args.parallel} threads/job")

    fun = arg_makers[args.function]
    it = iter_args(fun, rest, args.parallel, args.repeat, args.skip)  # pyright: ignore[reportArgumentType]
    if cli.dry_run:
        for cmd in it:
            print(" ".join(cmd))
    else:
        fs = [cli.thread_submit(subprocess.run, x, check=True) for x in it]
        cli.wait_raise(fs)
    print("End of " + __file__)


if __name__ == "__main__":
    main()
