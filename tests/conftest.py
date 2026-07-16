import os

# Some libfabric/OFI builds of meep's MPI backend mis-select a network
# provider on machines with an active VPN/utun interface, causing a clean
# simulation to abort with a nonzero exit code during MPI_Finalize at
# interpreter shutdown (results are unaffected — this is purely a finalize-
# time crash). Forcing the tcp provider avoids it; setdefault so a
# deliberately-set value is never overridden.
os.environ.setdefault("FI_PROVIDER", "tcp")

import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--runslow", action="store_true", default=False,
        help="run slow golden-value integration tests (real meep/meow sims)",
    )


def pytest_collection_modifyitems(config, items):
    if config.getoption("--runslow"):
        return
    skip_slow = pytest.mark.skip(reason="need --runslow option to run")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip_slow)
