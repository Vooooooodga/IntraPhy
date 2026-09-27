"""Task-local SQLite staging for exact, order-preserving match rows."""
from __future__ import annotations

from pathlib import Path
import pickle
import sqlite3


_BATCH_SIZE = 256
_READ_BATCH_SIZE = 128


def canonical_copy_pair_key(left, right):
    """Return a stable typed key for one unordered pair of gene copies."""
    return pickle.dumps(tuple(sorted((tuple(left), tuple(right)))), protocol=5)


class MatchStagingStore:
    """Retain full typed rows and prequalification facts in a caller-owned DB."""

    def __init__(self, path):
        self.path = Path(path)
        sidecars = (Path(str(self.path) + suffix) for suffix in ("-wal", "-shm", "-journal"))
        if self.path.exists() or self.path.is_symlink() or any(
            item.exists() or item.is_symlink() for item in sidecars
        ):
            raise FileExistsError(f"refusing to reuse match staging path: {self.path}")
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.connection = sqlite3.connect(self.path)
        try:
            self.connection.execute("PRAGMA journal_mode=DELETE")
            self.connection.execute("PRAGMA temp_store=FILE")
            self.connection.execute("PRAGMA cache_size=-8192")
            self.connection.execute("PRAGMA mmap_size=0")
            self.connection.execute(
                "CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL)"
            )
            self.connection.execute(
                "CREATE TABLE match_rows ("
                " ordinal INTEGER PRIMARY KEY, pair_key BLOB NOT NULL, payload BLOB NOT NULL,"
                " finalized INTEGER NOT NULL DEFAULT 0)"
            )
            self.connection.execute(
                "CREATE INDEX match_rows_pair_order ON match_rows(pair_key, ordinal)"
            )
            self.connection.execute(
                "INSERT INTO metadata(key,value) VALUES('schema','intraphy.match-stage/1')"
            )
            self.connection.execute(
                "INSERT INTO metadata(key,value) VALUES('status','building')"
            )
            self._pending = 0
            self.connection.commit()
        except BaseException as error:
            try:
                self.connection.rollback()
            finally:
                self.connection.close()
            try:
                failed = sqlite3.connect(self.path)
                try:
                    failed.execute(
                        "CREATE TABLE IF NOT EXISTS metadata "
                        "(key TEXT PRIMARY KEY, value TEXT NOT NULL)"
                    )
                    failed.execute(
                        "INSERT OR REPLACE INTO metadata(key,value) VALUES('schema',?)",
                        ("intraphy.match-stage/1",),
                    )
                    failed.execute(
                        "INSERT OR REPLACE INTO metadata(key,value) VALUES('status','incomplete')"
                    )
                    failed.execute(
                        "INSERT OR REPLACE INTO metadata(key,value) VALUES('detail',?)",
                        (f"{type(error).__name__}: {error}",),
                    )
                    failed.commit()
                finally:
                    failed.close()
            except BaseException:
                pass
            raise

    @staticmethod
    def _encode(row, projection_facts, overlap_pair):
        return pickle.dumps((row, projection_facts, overlap_pair), protocol=5)

    @staticmethod
    def _decode(payload):
        row, projection_facts, overlap_pair = pickle.loads(payload)
        return row, projection_facts, overlap_pair

    def add(self, ordinal, pair_key, row, projection_facts=(), overlap_pair=None):
        payload = self._encode(row, tuple(projection_facts), overlap_pair)
        self.connection.execute(
            "INSERT INTO match_rows(ordinal,pair_key,payload,finalized) VALUES(?,?,?,0)",
            (int(ordinal), sqlite3.Binary(pair_key), sqlite3.Binary(payload)),
        )
        self._pending += 1
        if self._pending >= _BATCH_SIZE:
            self.connection.commit()
            self._pending = 0

    def pair_keys(self):
        self.connection.commit()
        keys = self.connection.execute(
            "SELECT pair_key FROM match_rows GROUP BY pair_key ORDER BY MIN(ordinal)"
        ).fetchall()
        for (key,) in keys:
            yield bytes(key)

    def pair_count(self):
        return int(self.connection.execute(
            "SELECT COUNT(DISTINCT pair_key) FROM match_rows"
        ).fetchone()[0])

    def row_count(self):
        return int(self.connection.execute("SELECT COUNT(*) FROM match_rows").fetchone()[0])

    def pair_rows(self, pair_key):
        cursor = self.connection.execute(
            "SELECT ordinal,payload FROM match_rows WHERE pair_key=? ORDER BY ordinal",
            (sqlite3.Binary(pair_key),),
        )
        return [(ordinal, *self._decode(payload)) for ordinal, payload in cursor]

    def replace_pair_rows(self, pair_key, rows):
        for ordinal, row, projection_facts, overlap_pair in rows:
            payload = self._encode(row, projection_facts, overlap_pair)
            self.connection.execute(
                "UPDATE match_rows SET payload=?, finalized=1 WHERE ordinal=? AND pair_key=?",
                (sqlite3.Binary(payload), int(ordinal), sqlite3.Binary(pair_key)),
            )
        self.connection.commit()

    def replace_row(self, ordinal, row, projection_facts=(), overlap_pair=None):
        payload = self._encode(row, tuple(projection_facts), overlap_pair)
        self.connection.execute(
            "UPDATE match_rows SET payload=?, finalized=1 WHERE ordinal=?",
            (sqlite3.Binary(payload), int(ordinal)),
        )
        self._pending += 1
        if self._pending >= _BATCH_SIZE:
            self.connection.commit()
            self._pending = 0

    def iter_rows(self):
        self.connection.commit()
        last_ordinal = -1
        while True:
            rows = self.connection.execute(
                "SELECT ordinal,payload FROM match_rows WHERE ordinal>? "
                "ORDER BY ordinal LIMIT ?",
                (last_ordinal, _READ_BATCH_SIZE),
            ).fetchall()
            if not rows:
                break
            for ordinal, payload in rows:
                last_ordinal = ordinal
                yield (ordinal, *self._decode(payload))

    def set_status(self, status, detail=""):
        if status not in {"building", "resolving", "aggregating", "complete", "incomplete"}:
            raise ValueError(f"invalid match staging status: {status}")
        self.connection.execute(
            "INSERT OR REPLACE INTO metadata(key,value) VALUES('status',?)", (status,)
        )
        if detail:
            self.connection.execute(
                "INSERT OR REPLACE INTO metadata(key,value) VALUES('detail',?)", (detail,)
            )
        self.connection.commit()

    def close(self):
        self.connection.rollback()
        self.connection.close()
