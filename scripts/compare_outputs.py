"""Compare two scaffold output trees (e.g. output_omnium vs output_remake) file by file.

Each file present in either tree is classified as:
    identical  - same values (and same structure)
    close      - same structure, values within tolerance: |a - b| <= atol + rtol * |a|,
                 with atol = atol_frac * max|a| for that variable
    DIFFERENT  - structure differs, or values outside tolerance
    MISSING_A / MISSING_B - only present in one tree
    ERROR      - could not be read/compared
    present    - figures (png/pdf): only presence is checked

Usage:
    python compare_outputs.py <tree_a> <tree_b> [--rtol 1e-3] [--atol-frac 1e-6]
        [--include REGEX] [--exclude REGEX] [-j NPROC] [--report report.csv]

Runs in either pixi env (omnium: py3.7, remake: py3.12).
"""
import argparse
import csv
import os
import pickle
import re
import sys
from collections import OrderedDict
from multiprocessing import Pool

import numpy as np

SKIP_RE = re.compile(r'(\.done$|\.log$|(^|/)metadata/|\.name$|^settings\.json$|\.tgz$)')
FIGURE_EXTS = {'.png', '.pdf', '.svg'}
NUM_RE = re.compile(r'^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$')
IGNORE_ATTR_RE = re.compile(r'^(omnium_(vn|git_\w+)|\w+_vn|\w+_git_\w+|history)$')

# Older cloud_tracking pickles used this module path (renamed in cloud_tracking 8b240fa).
try:
    import cloud_tracking.tracking
    sys.modules.setdefault('cloud_tracking.cloud_tracking', cloud_tracking.tracking)
except ImportError:
    pass


class Result:
    def __init__(self, status, detail=''):
        self.status = status
        self.detail = detail


def norm_name(name):
    # iris/netCDF versions differ in sanitising variable names (w1.0 vs w1_0).
    return re.sub(r'[.\-]', '_', name)


def compare_arrays(a, b, rtol, atol_frac):
    """Return (status, detail) for two numeric arrays."""
    a = np.ma.asarray(a)
    b = np.ma.asarray(b)
    if a.shape != b.shape:
        return 'DIFFERENT', 'shape {} vs {}'.format(a.shape, b.shape)
    ma, mb = np.ma.getmaskarray(a), np.ma.getmaskarray(b)
    if not np.array_equal(ma, mb):
        return 'DIFFERENT', 'mask differs at {} points'.format(int((ma != mb).sum()))
    if a.dtype.kind not in 'biufc' or b.dtype.kind not in 'biufc':
        same = np.array_equal(np.ma.getdata(a), np.ma.getdata(b))
        return ('identical', '') if same else ('DIFFERENT', 'non-numeric values differ')
    da = np.ma.getdata(a)[~ma].astype(float)
    db = np.ma.getdata(b)[~mb].astype(float)
    nan_a, nan_b = np.isnan(da), np.isnan(db)
    if not np.array_equal(nan_a, nan_b):
        return 'DIFFERENT', 'NaN pattern differs at {} points'.format(int((nan_a != nan_b).sum()))
    da, db = da[~nan_a], db[~nan_b]
    if np.array_equal(da, db):
        return 'identical', ''
    diff = np.abs(da - db)
    scale = np.max(np.abs(da)) if da.size else 0.
    atol = atol_frac * scale
    bad = diff > atol + rtol * np.abs(da)
    with np.errstate(divide='ignore', invalid='ignore'):
        rel = np.where(da != 0, diff / np.abs(da), np.where(diff > 0, np.inf, 0.))
    detail = 'max_abs {:.3g} max_rel {:.3g} n_diff {}/{}'.format(
        diff.max(), np.nanmax(rel), int((diff > 0).sum()), diff.size)
    if bad.any():
        return 'DIFFERENT', detail + ' n_out_of_tol {}'.format(int(bad.sum()))
    return 'close', detail


def compare_labels(a, b):
    """Labelled-cloud fields: labels may be numbered differently but must partition alike."""
    a = np.ma.getdata(np.ma.asarray(a)).astype(np.int64)
    b = np.ma.getdata(np.ma.asarray(b)).astype(np.int64)
    if a.shape != b.shape:
        return 'DIFFERENT', 'shape {} vs {}'.format(a.shape, b.shape)
    if np.array_equal(a, b):
        return 'identical', ''
    if not np.array_equal(a > 0, b > 0):
        return 'DIFFERENT', 'cloud mask differs at {} points'.format(int(((a > 0) != (b > 0)).sum()))
    pairs = np.unique(np.stack([a[a > 0], b[b > 0]]), axis=1)
    if len(np.unique(pairs[0])) == pairs.shape[1] == len(np.unique(pairs[1])):
        return 'close', 'labels renumbered (same partition)'
    return 'DIFFERENT', 'cloud partition differs'


def worst(statuses):
    order = ['identical', 'close', 'present', 'DIFFERENT', 'ERROR']
    return max(statuses, key=order.index) if statuses else 'identical'


def compare_nc(fa, fb, rtol, atol_frac):
    import netCDF4
    details = []
    statuses = []
    with netCDF4.Dataset(fa) as a, netCDF4.Dataset(fb) as b:
        va = OrderedDict((norm_name(n), n) for n in a.variables)
        vb = OrderedDict((norm_name(n), n) for n in b.variables)
        only = sorted(set(va) ^ set(vb))
        if only:
            statuses.append('DIFFERENT')
            details.append('vars only in one: {}'.format(only[:10]))
        for k in va:
            if k not in vb:
                continue
            xa, xb = a.variables[va[k]], b.variables[vb[k]]
            xa.set_auto_mask(True)
            xb.set_auto_mask(True)
            if k.startswith('labelled_clouds'):
                st, d = compare_labels(xa[:], xb[:])
            else:
                st, d = compare_arrays(xa[:], xb[:], rtol, atol_frac)
            attrs_a = {n: str(xa.getncattr(n)) for n in xa.ncattrs() if not IGNORE_ATTR_RE.match(n)}
            attrs_b = {n: str(xb.getncattr(n)) for n in xb.ncattrs() if not IGNORE_ATTR_RE.match(n)}
            # Variable-name references (e.g. coordinates) may be sanitised differently.
            attrs_a = {n: norm_name(v) for n, v in attrs_a.items()}
            attrs_b = {n: norm_name(v) for n, v in attrs_b.items()}
            if attrs_a != attrs_b:
                diff_keys = sorted(n for n in set(attrs_a) | set(attrs_b)
                                   if attrs_a.get(n) != attrs_b.get(n))
                d = (d + '; ' if d else '') + 'attrs differ: {}'.format(diff_keys)
                st = worst([st, 'close'])
            statuses.append(st)
            if st != 'identical':
                details.append('{}: {}'.format(k, d))
    return Result(worst(statuses), ' | '.join(details))


def compare_frames(da, db, rtol, atol_frac):
    if list(da.columns) != list(db.columns) or da.shape != db.shape:
        return 'DIFFERENT', 'columns/shape {} {} vs {} {}'.format(
            list(da.columns), da.shape, list(db.columns), db.shape)
    statuses, details = [], []
    for col in da.columns:
        ca, cb = da[col].values, db[col].values
        if ca.dtype.kind in 'biufc' and cb.dtype.kind in 'biufc':
            st, d = compare_arrays(ca, cb, rtol, atol_frac)
        else:
            same = [str(x) for x in ca] == [str(x) for x in cb]
            st, d = ('identical', '') if same else ('DIFFERENT', 'values differ')
        statuses.append(st)
        if st != 'identical':
            details.append('{}: {}'.format(col, d))
    return worst(statuses), ' | '.join(details)


def compare_hdf(fa, fb, rtol, atol_frac):
    import pandas as pd
    with pd.HDFStore(fa, 'r') as sa, pd.HDFStore(fb, 'r') as sb:
        ka, kb = sorted(sa.keys()), sorted(sb.keys())
        if ka != kb:
            return Result('DIFFERENT', 'keys {} vs {}'.format(ka, kb))
        statuses, details = [], []
        for k in ka:
            st, d = compare_frames(sa[k], sb[k], rtol, atol_frac)
            statuses.append(st)
            if d:
                details.append('{}: {}'.format(k, d))
    return Result(worst(statuses), ' | '.join(details))


def compare_csv(fa, fb, rtol, atol_frac):
    import pandas as pd
    st, d = compare_frames(pd.read_csv(fa), pd.read_csv(fb), rtol, atol_frac)
    return Result(st, d)


def compare_text(fa, fb, rtol, atol_frac):
    with open(fa) as f:
        ta = f.read()
    with open(fb) as f:
        tb = f.read()
    if ta == tb:
        return Result('identical')
    toks_a = re.split(r'([\s,&\\]+)', ta)
    toks_b = re.split(r'([\s,&\\]+)', tb)
    if len(toks_a) != len(toks_b):
        return Result('DIFFERENT', 'different number of tokens')
    na, nb = [], []
    for x, y in zip(toks_a, toks_b):
        if NUM_RE.match(x) and NUM_RE.match(y):
            na.append(float(x))
            nb.append(float(y))
        elif x != y:
            return Result('DIFFERENT', 'text differs: {!r} vs {!r}'.format(x[:40], y[:40]))
    st, d = compare_arrays(np.array(na), np.array(nb), rtol, atol_frac)
    return Result(st, d)


def tracker_summary(tracker):
    """Comparable summary of a cloud_tracking Tracker (the object graph is cyclic)."""
    groups = sorted((int(g.is_linear), len(g.clds),
                     tuple(sorted(c.lifetime for c in g.end_clouds)))
                    for g in tracker.groups)
    return {'n_clds': len(tracker.all_clds), 'n_groups': len(tracker.groups), 'groups': groups}


def compare_objects(a, b, rtol, atol_frac, path='obj'):
    """Recursively compare unpickled objects. Returns (status, details list)."""
    if type(a).__name__ == 'Tracker' and type(b).__name__ == 'Tracker':
        sa, sb = tracker_summary(a), tracker_summary(b)
        return ('identical', []) if sa == sb else (
            'DIFFERENT', ['{}: tracker summary differs (n_clds {} vs {}, n_groups {} vs {})'.format(
                path, sa['n_clds'], sb['n_clds'], sa['n_groups'], sb['n_groups'])])
    if isinstance(a, dict) and isinstance(b, dict):
        if set(map(str, a)) != set(map(str, b)):
            return 'DIFFERENT', ['{}: keys differ'.format(path)]
        statuses, details = [], []
        bmap = {str(k): v for k, v in b.items()}
        for k, v in a.items():
            st, d = compare_objects(v, bmap[str(k)], rtol, atol_frac, '{}[{!r}]'.format(path, k))
            statuses.append(st)
            details.extend(d)
        return worst(statuses), details
    if isinstance(a, (list, tuple)) and isinstance(b, (list, tuple)):
        if len(a) != len(b):
            return 'DIFFERENT', ['{}: length {} vs {}'.format(path, len(a), len(b))]
        if all(isinstance(x, (int, float, np.number)) for x in list(a) + list(b)):
            st, d = compare_arrays(np.array(a), np.array(b), rtol, atol_frac)
            return st, ['{}: {}'.format(path, d)] if d else []
        statuses, details = [], []
        for i, (x, y) in enumerate(zip(a, b)):
            st, d = compare_objects(x, y, rtol, atol_frac, '{}[{}]'.format(path, i))
            statuses.append(st)
            details.extend(d)
        return worst(statuses), details
    if isinstance(a, (np.ndarray, np.number, int, float)) and \
            isinstance(b, (np.ndarray, np.number, int, float)):
        st, d = compare_arrays(np.atleast_1d(a), np.atleast_1d(b), rtol, atol_frac)
        return st, ['{}: {}'.format(path, d)] if d else []
    if hasattr(a, '__dict__') and hasattr(b, '__dict__') and type(a).__name__ == type(b).__name__:
        return compare_objects(vars(a), vars(b), rtol, atol_frac, path)
    return ('identical', []) if a == b else ('DIFFERENT', ['{}: {!r} vs {!r}'.format(
        path, str(a)[:40], str(b)[:40])])


def compare_pickle(fa, fb, rtol, atol_frac):
    sys.setrecursionlimit(100000)
    with open(fa, 'rb') as f:
        a = pickle.load(f)
    with open(fb, 'rb') as f:
        b = pickle.load(f)
    st, details = compare_objects(a, b, rtol, atol_frac)
    return Result(st, ' | '.join(details[:10]))


COMPARERS = {
    '.nc': compare_nc,
    '.hdf': compare_hdf,
    '.csv': compare_csv,
    '.txt': compare_text,
    '.tex': compare_text,
    '.pkl': compare_pickle,
    '.np': compare_pickle,
    '.dummy': compare_text,
}


def compare_file(args):
    rel, tree_a, tree_b, rtol, atol_frac = args
    fa, fb = os.path.join(tree_a, rel), os.path.join(tree_b, rel)
    if not os.path.exists(fa):
        return rel, Result('MISSING_A')
    if not os.path.exists(fb):
        return rel, Result('MISSING_B')
    ext = os.path.splitext(rel)[1]
    if ext in FIGURE_EXTS:
        return rel, Result('present')
    comparer = COMPARERS.get(ext, compare_text)
    with open(fa, 'rb') as f:
        if f.read(1) == b'\x80':
            # A pickle whatever the extension (e.g. precip_plot's *_cached_hov_*.nc are
            # written with ndarray.dump).
            comparer = compare_pickle
    try:
        return rel, comparer(fa, fb, rtol, atol_frac)
    except Exception as e:
        return rel, Result('ERROR', '{}: {}'.format(type(e).__name__, str(e)[:200]))


def list_files(tree, include, exclude):
    rels = []
    for root, dirs, files in os.walk(tree):
        dirs.sort()
        for fn in files:
            rel = os.path.relpath(os.path.join(root, fn), tree)
            if SKIP_RE.search(rel):
                continue
            if include and not re.search(include, rel):
                continue
            if exclude and re.search(exclude, rel):
                continue
            rels.append(rel)
    return rels


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('tree_a')
    parser.add_argument('tree_b')
    parser.add_argument('--rtol', type=float, default=1e-3)
    parser.add_argument('--atol-frac', type=float, default=1e-6)
    parser.add_argument('--include')
    parser.add_argument('--exclude')
    parser.add_argument('-j', '--nproc', type=int, default=1)
    parser.add_argument('--report', help='write a CSV with one row per file')
    args = parser.parse_args()

    rels = sorted(set(list_files(args.tree_a, args.include, args.exclude)) |
                  set(list_files(args.tree_b, args.include, args.exclude)))
    tasks = [(rel, args.tree_a, args.tree_b, args.rtol, args.atol_frac) for rel in rels]
    if args.nproc > 1:
        with Pool(args.nproc) as pool:
            results = pool.map(compare_file, tasks, chunksize=1)
    else:
        results = [compare_file(t) for t in tasks]

    counts = OrderedDict()
    for rel, res in results:
        counts[res.status] = counts.get(res.status, 0) + 1
        if res.status not in ('identical', 'present'):
            print('{:10s} {}  {}'.format(res.status, rel, res.detail))
    print('\nsummary ({} files, rtol={}, atol_frac={}): {}'.format(
        len(results), args.rtol, args.atol_frac,
        ', '.join('{} {}'.format(k, v) for k, v in counts.items())))

    if args.report:
        with open(args.report, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['path', 'status', 'detail'])
            for rel, res in results:
                writer.writerow([rel, res.status, res.detail])

    bad = {'DIFFERENT', 'MISSING_A', 'MISSING_B', 'ERROR'}
    return 1 if bad & set(counts) else 0


if __name__ == '__main__':
    sys.exit(main())
