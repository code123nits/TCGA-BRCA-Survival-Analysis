"""Verify local processed inputs against archived SHA-256 hashes."""
import argparse, csv, hashlib
from pathlib import Path

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data', type=Path, required=True, help='Directory holding three processed input files')
    ap.add_argument('--manifest', type=Path, default=Path(__file__).resolve().parents[1] / 'documentation/input_sha256_manifest.csv')
    args = ap.parse_args()
    ok = True
    with args.manifest.open(newline='') as f:
        for row in csv.DictReader(f):
            path = args.data / row['filename']
            if not path.is_file():
                print(f'MISSING: {path}'); ok = False; continue
            h = hashlib.sha256()
            with path.open('rb') as source:
                for chunk in iter(lambda: source.read(4 * 1024 * 1024), b''):
                    h.update(chunk)
            match = h.hexdigest() == row['sha256'] and path.stat().st_size == int(row['size_bytes'])
            print(('PASS' if match else 'FAIL') + ': ' + row['filename'])
            ok &= match
    if not ok: raise SystemExit(1)
    print('ALL INPUT HASH CHECKS PASSED')

if __name__ == '__main__': main()
