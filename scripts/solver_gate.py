"""Required complete-solve comparison; never skip missing or unconverged cases."""
import argparse
import json
import math
import re
from pathlib import Path
import subprocess

from flexfoil import _rustfoil as native

ROOT = Path(__file__).resolve().parents[1]
# XFOIL polar print half-units: CL/CM/trips 5e-5, CD 5e-6; both solvers stop at 1e-4 RMS.
# These are bounded regression tolerances for this matrix, not physical-accuracy claims.
TOLERANCES = dict(cl=2e-4, cd=1e-5, cd_pressure=1e-5, cd_friction=1e-5, cm=1e-4, x_tr_upper=2e-4, x_tr_lower=2e-4)


def require_close(actual, expected, tolerance, label):
    if not math.isfinite(actual) or not math.isfinite(expected) or abs(actual - expected) > tolerance:
        raise AssertionError(f'{label}: {actual} vs {expected}, tolerance {tolerance}')


def run_reference(binary, geometry, alpha, reynolds, mach, ncrit, trips, directory):
    (directory / 'foil.dat').write_bytes(geometry.read_bytes())
    # PCOP preserves identical points. In-memory PACC/PWRT avoids legacy EOF append handling.
    commands = (f'PLOP\nG\n\nLOAD foil.dat\nPCOP\nOPER\nVISC {reynolds:g}\nMACH {mach:g}\n'
                f'ITER 100\nVPAR\nN {ncrit:g}\nXTR {trips[0]:g} {trips[1]:g}\n\n'
                f'PACC\n\n\nALFA {alpha:g}\nCPWR cp.txt\nDUMP bl.txt\nPWRT\npolar.txt\n\nQUIT\n')
    (directory / 'commands.in').write_text(commands)
    process = subprocess.run([str(binary)], input=commands, text=True, capture_output=True,
                             cwd=directory, timeout=60)
    (directory / 'stdout.log').write_text(process.stdout)
    (directory / 'stderr.log').write_text(process.stderr)
    if process.returncode:
        raise AssertionError(f'XFOIL exit {process.returncode}: {process.stderr}')
    lines = (directory / 'polar.txt').read_text().splitlines()
    text = '\n'.join(lines)
    settings = re.search(r'Mach =\s*([\d.]+)\s+Re =\s*([\d.]+)\s*e\s*([+-]?\d+)\s+Ncrit =\s*([\d.]+)\s+([\d.]+)', text)
    trip_settings = re.search(r'xtrf =\s*([\d.]+) \(top\)\s*([\d.]+)', text)
    if settings is None or trip_settings is None or '1 1 Reynolds number fixed' not in text:
        raise AssertionError('Reference operating-condition header missing')
    m, r, exponent, n1, n2 = map(float, settings.groups())
    for actual, expected in [(m, mach), (r * 10**exponent, reynolds), (n1, ncrit), (n2, ncrit), *zip(map(float, trip_settings.groups()), trips)]:
        require_close(actual, expected, 1e-8, 'reference operating condition')
    separator = next(i for i, line in enumerate(lines) if line.strip().startswith('------'))
    values = [list(map(float, line.split())) for line in lines[separator + 1:] if line.strip()]
    if len(values) != 1 or len(values[0]) != 9:
        raise AssertionError(f'XFOIL must write exactly one converged polar row: {values}')
    row = values[0]
    require_close(row[0], alpha, 1e-8, 'reference alpha')
    result = dict(zip(['alpha', 'cl', 'cd', 'cdp', 'cm', 'x_tr_upper', 'x_tr_lower'], row))
    # XFOIL polar CDp is the surface-pressure integral, not its reported CD-CDf component.
    components = re.findall(r'CDf =\s*([+\-\d.eE]+)\s+CDp =\s*([+\-\d.eE]+)', process.stdout)
    if not components:
        raise AssertionError('Reference drag components missing')
    result['cd_friction'], result['cd_pressure'] = map(float, components[-1])
    distribution = [list(map(float, line.split())) for line in
                    (directory / 'bl.txt').read_text().splitlines()[1:] if line.strip()]
    points = [list(map(float, line.split())) for line in geometry.read_text().splitlines()[1:]]
    if len(distribution) < len(points):
        raise AssertionError('Reference surface stations missing')
    for point, station in zip(points, distribution):
        for actual, expected in zip(station[1:3], point):
            require_close(actual, expected, 5.1e-6, 'reference panel coordinate')
    return result, distribution[:len(points)]


def compare_distributions(bl, distribution, mach):
    if not bl.get('success') or not bl.get('converged'):
        raise AssertionError('Boundary-layer replay did not converge')
    le = min(range(len(distribution)), key=lambda i: distribution[i][1])
    checked = 0
    for side, reference in [('upper', distribution[:le + 1]), ('lower', distribution[le + 1:])]:
        xs = bl['x_' + side]
        for key in ('ue', 'theta', 'delta_star', 'h', 'cf'):
            if len(bl[key + '_' + side]) != len(xs):
                raise AssertionError('Boundary-layer array lengths differ')
        for i, x in enumerate(xs):
            # The synthetic stagnation row and cross-LE nose mapping are a separate station test.
            if x < .05:
                continue
            row = min(reference, key=lambda point: abs(point[1] - x))
            require_close(x, row[1], 5.1e-6, 'BL panel location')
            ui = bl['ue_' + side][i]
            beta = math.sqrt(1 - mach * mach)
            tk = mach * mach / (1 + beta)**2
            ue = ui * (1 - tk) / (1 - tk * ui * ui)
            local_msq = ue * ue * mach * mach / (1 + .2 * mach * mach * (1 - ue * ue))
            density = ((1 + .2 * mach * mach) / (1 + .2 * local_msq))**2.5
            h = bl['h_' + side][i]
            hk = (h - .29 * local_msq) / (1 + .113 * local_msq)
            # xoper.f BLDUMP reports local Cf * rho * Ue**2 and kinematic Hk.
            values = [(ue, abs(row[3]), 2e-4),
                      (bl['delta_star_' + side][i], row[4], 5e-6),
                      (bl['theta_' + side][i], row[5], 2e-6),
                      (bl['cf_' + side][i] * density * ue**2, row[6], 2e-5),
                      (hk, row[7], .005)]
            for key, (actual, expected, tolerance) in zip(('ue', 'dstar', 'theta', 'cf', 'hk'), values):
                require_close(actual, expected, tolerance, f'{side} {key} at x={x}')
            checked += 1
    if checked < 100:
        raise AssertionError(f'Too few matched BL stations: {checked}')
    return checked


def gate(binary, output, mach_values):
    output.mkdir(parents=True, exist_ok=True)
    cases = []
    failures = []
    for foil in ('0012', '2412', '4412'):
        geometry = ROOT / f'testdata/solver-reference/naca{foil}.dat'
        points = [list(map(float, line.split())) for line in geometry.read_text().splitlines()[1:]]
        flat = [value for point in points for value in point]
        for reynolds in (1e5, 1e6):
            for mach in mach_values:
                serial = []
                for alpha, ncrit, trips in ((0., 9., (1., 1.)), (4., 9., (1., 1.)), (8., 9., (1., 1.)), (4., 5., (.3, .4))):
                    name = f'naca{foil}-re{reynolds:g}-a{alpha:g}-m{mach:g}-n{ncrit:g}-tr{trips[0]:g}-{trips[1]:g}'
                    directory = output / name
                    directory.mkdir(exist_ok=True)
                    for file in ('polar.txt', 'cp.txt', 'bl.txt'):
                        (directory / file).unlink(missing_ok=True)
                    reference, distribution = run_reference(binary, geometry, alpha, reynolds, mach, ncrit, trips, directory)
                    rust = native.analyze_faithful(flat, alpha, reynolds, mach, ncrit, 100, 1, *trips)
                    serial.append(rust)
                    bl = native.get_bl_distribution(flat, alpha, reynolds, mach, ncrit, 100, 1, *trips)
                    case = dict(name=name, coords=flat, alpha=alpha, reynolds=reynolds, mach=mach,
                                ncrit=ncrit, trips=trips, reference=reference, rust=rust, bl=bl,
                                cp_reference=[list(map(float, line.split())) for line in (directory / 'cp.txt').read_text().splitlines()[1:] if line.strip()])
                    cases.append(case)
                    try:
                        if not rust.get('success') or not rust.get('converged'):
                            raise AssertionError(f'Rust did not converge: {rust}')
                        require_close(rust['residual'], 0., 1e-4, 'converged RMS')
                        if not 1 <= rust['iterations'] <= 100:
                            raise AssertionError('Invalid converged iteration count')
                        case['matched_bl_stations'] = compare_distributions(bl, distribution, mach)
                        for key, tolerance in TOLERANCES.items():
                            require_close(rust[key], reference[key], tolerance, key)
                    except AssertionError as error:
                        failures.append(f'{name}: {error}')
                    print(name, 'converged', rust.get('converged'), flush=True)
                batch = native.analyze_faithful_batch(flat, [0., 4., 8.], reynolds, mach, 9., 100)
                if len(batch) != 3 or len(serial) != 4:
                    raise AssertionError('Batch dropped cases')
                for first, second in zip(serial[:3], batch):
                    if first.get('converged') != second.get('converged') or first.get('success') != second.get('success'):
                        raise AssertionError('Serial/batch status mismatch')
                    for key in TOLERANCES:
                        require_close(first[key], second[key], 1e-12, 'serial/batch ' + key)
    if len(cases) != 24 * len(mach_values):
        raise AssertionError('Required cases missing')
    (output / 'report.json').write_text(json.dumps(dict(cases=cases, failures=failures), indent=2) + '\n')
    if failures:
        raise AssertionError('\n'.join(failures))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--xfoil', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--mach-zero-only', action='store_true', help='diagnostic subset; never used by CI')
    args = parser.parse_args()
    gate(args.xfoil.resolve(strict=True), args.output.resolve(), [0.] if args.mach_zero_only else [0., .25])
