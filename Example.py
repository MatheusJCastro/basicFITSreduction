"""Backward-compatible entry point for the ReduceRaw package."""

from reduce_raw import ReduceRaw


def main() -> None:
    """Run the example configuration defined for this project."""
    reduction = ReduceRaw()
    reduction.dir_name = "2026-10-06/group2"
    reduction.obj_name = "NGC7635"
    reduction.combine_filters = ["R", "G", "B"]
    reduction.combine = True
    reduction.reduce_bias = True
    reduction.reduce_flat = True
    reduction.run()


if __name__ == "__main__":
    main()
