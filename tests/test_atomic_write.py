"""Can a concurrent reader ever see a half-written file? (dc_validation.md item 6)

Every SLURM array sweep in examples/ ends by calling `aggregate()`, which reads
*every* task's per-task output -- including ones a peer task is writing at that
instant. Those per-task files were written with a plain `np.savez`, which is
not atomic: a reader that opens the path between the writer's open() and its
last flush gets a truncated zip.

MEASURED. Job 11517245 (the calibrated stub sweep) lost **7 of its 9 tasks** to

    zipfile.BadZipFile: File is not a zip file

raised from `np.load` inside `aggregate()`. The run's result survived only by
luck: two tasks happened to aggregate after every write had finished. The same
code had passed in job 11500412 purely because those tasks finished 12 minutes
apart instead of 6 -- a latent race, not a fixed one.

The fix is `beef.atomic_replace_write`: write to a temp path, then os.replace,
which is atomic on POSIX so a reader sees the old file or the new one and never
a partial one. These tests guard the helper, not the sweeps.

TWO TRAPS ENCODED HERE, both of which have already cost a run
-------------------------------------------------------------
1. **`np.savez` given a path without a `.npz` suffix appends one.** Handed the
   temp name it would write `<tmp>.npz`, leave `<tmp>` absent, and the rename
   would fail *after* the simulation work was done -- which is exactly how all
   four tasks of an earlier job died (dc_validation.md Done item 6). Callers
   must pass an open file object. The helper now raises with that diagnosis
   rather than letting a bare FileNotFoundError surface.

2. **The temp name must carry the hostname, not just the PID.** PIDs are
   allocated per host, so a PID-only temp name is not unique on a shared
   filesystem spanning nodes (dc_validation.md item 4). The seven example-local
   copies of this helper all used the weaker PID-only form; the library one
   does not, and `test_tmp_name_is_host_qualified` keeps it that way.
"""

import os
import socket
import threading
import time

import numpy as np
import pytest

from beef import atomic_replace_write


def _writer(payload):
    """The documented-correct write_fn: file object, closed via `with`.

    Not `lambda t: np.savez(open(t, "wb"), ...)` — zipfile does not close a
    file object it was handed, so that form leaves the flush to refcounting.
    It happens to work on CPython and is exactly the kind of thing that stops
    working somewhere else, so the tests use the form the sweeps use."""
    def write_fn(tmp):
        with open(tmp, "wb") as fh:
            np.savez(fh, a=payload)
    return write_fn


class TestTheFileIsCompleteWhenItAppears:
    def test_target_holds_what_was_written(self, tmp_path):
        path = str(tmp_path / "out.npz")
        payload = np.arange(1000, dtype=float)
        atomic_replace_write(path, _writer(payload))
        assert np.array_equal(np.load(path)["a"], payload)

    def test_no_temp_residue_left_behind(self, tmp_path):
        path = str(tmp_path / "out.npz")
        atomic_replace_write(path, _writer(np.zeros(3)))
        leftovers = [p for p in os.listdir(tmp_path) if ".tmp." in p]
        assert leftovers == [], f"temp files not cleaned up: {leftovers}"

    def test_overwriting_an_existing_file_is_also_atomic(self, tmp_path):
        # The sweeps re-run indices, so the target usually already exists.
        path = str(tmp_path / "out.npz")
        atomic_replace_write(path, _writer(np.ones(3)))
        atomic_replace_write(path, _writer(np.full(3, 7.0)))
        assert np.array_equal(np.load(path)["a"], np.full(3, 7.0))


class TestTheNpzSuffixTrap:
    """np.savez(path_without_npz) writes somewhere else. Twice now, that has
    destroyed a job's worth of completed simulation."""

    def test_string_form_is_rejected_with_a_diagnosis(self, tmp_path):
        path = str(tmp_path / "out.npz")
        with pytest.raises(RuntimeError, match="np.savez"):
            atomic_replace_write(path, lambda t: np.savez(t, a=np.zeros(3)))

    def test_the_target_is_not_created_by_the_failed_attempt(self, tmp_path):
        # Better to have no file than a file the aggregate would misread.
        path = str(tmp_path / "out.npz")
        with pytest.raises(RuntimeError):
            atomic_replace_write(path, lambda t: np.savez(t, a=np.zeros(3)))
        assert not os.path.exists(path)

    def test_file_object_form_is_the_documented_fix(self, tmp_path):
        path = str(tmp_path / "out.npz")
        atomic_replace_write(path, _writer(np.zeros(3)))
        assert os.path.exists(path)


class TestTheTempNameIsHostQualified:
    """dc_validation.md item 4 -- a PID-only temp name is not unique across
    nodes, and these jobs run 48 ranks over 2 of them."""

    def test_tmp_name_carries_hostname_and_pid(self, tmp_path):
        seen = []
        atomic_replace_write(str(tmp_path / "out.bin"),
                             lambda t: (seen.append(t), open(t, "wb").write(b"x")))
        assert len(seen) == 1
        name = os.path.basename(seen[0])
        assert socket.gethostname() in name, f"hostname missing from temp name {name!r}"
        assert str(os.getpid()) in name, f"pid missing from temp name {name!r}"


class TestAConcurrentReaderNeverSeesAPartialFile:
    """The actual failure mode, reproduced: a reader polling the path while a
    writer works must get either nothing, or a complete loadable file."""

    def test_reader_never_observes_a_truncated_zip(self, tmp_path):
        path = str(tmp_path / "out.npz")
        # Big enough that the write is not instantaneous -- a small payload
        # would pass this test even with a plain savez.
        payload = np.random.default_rng(0).random((200, 500))

        errors, reads = [], []
        stop = threading.Event()

        def poll():
            while not stop.is_set():
                if os.path.exists(path):
                    try:
                        reads.append(np.load(path)["a"].shape)
                    except Exception as exc:      # BadZipFile is the one that bit us
                        errors.append(exc)
                time.sleep(0.0005)

        reader = threading.Thread(target=poll, daemon=True)
        reader.start()
        try:
            for _ in range(20):
                atomic_replace_write(path, _writer(payload))
        finally:
            stop.set()
            reader.join(timeout=5)

        assert not errors, f"reader saw a partial file: {errors[:3]}"
        assert reads, "reader never managed to open the file — test proves nothing"
        assert all(s == payload.shape for s in reads)

    def test_the_same_probe_catches_a_plain_write(self, tmp_path):
        """Control: the polling probe above is sensitive enough to detect the
        bug it is guarding against. Without this, a probe that never races
        would pass the previous test vacuously."""
        path = str(tmp_path / "plain.npz")
        payload = np.random.default_rng(1).random((200, 500))

        errors = []
        stop = threading.Event()

        def poll():
            while not stop.is_set():
                if os.path.exists(path):
                    try:
                        np.load(path)["a"]
                    except Exception as exc:
                        errors.append(type(exc).__name__)
                time.sleep(0.0005)

        reader = threading.Thread(target=poll, daemon=True)
        reader.start()
        try:
            for _ in range(40):
                np.savez(path, a=payload)          # the unsafe pattern
        finally:
            stop.set()
            reader.join(timeout=5)

        # Timing-dependent by nature, so this is not asserted as a hard
        # requirement -- it is recorded so a future reader can see the probe
        # does race in practice. If this ever stops catching anything the
        # test above has quietly become vacuous.
        if not errors:
            pytest.skip("probe did not hit the race on this machine this run")
