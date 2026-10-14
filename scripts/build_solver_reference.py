"""Build the pinned XFOIL 6.99 oracle; missing tools/downloads are errors."""
import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile

SHA256 = '5c0250643f52ce0e75d7338ae2504ce7907f2d49a30f921826717b8ac12ebe40'
URL = 'https://web.mit.edu/drela/Public/web/xfoil/xfoil6.99.tgz'


def build(output):
    output.mkdir(parents=True, exist_ok=True)
    archive = output / 'xfoil6.99.tgz'
    if not archive.exists():
        subprocess.run(['curl', '--fail', '--location', '--proto', '=https',
                        '--max-time', '120', URL, '--output', str(archive)], check=True)
    if hashlib.sha256(archive.read_bytes()).hexdigest() != SHA256:
        raise RuntimeError('XFOIL source checksum mismatch')
    with tarfile.open(archive) as source:
        source.extractall(output, filter='data')
    compiler = os.environ.get('FC', 'gfortran')
    flags = '-O2 -fdefault-real-8 -fdefault-double-8 -fallow-argument-mismatch'
    cflags = subprocess.check_output(['pkg-config', '--cflags', 'x11'], text=True).strip()
    libraries = subprocess.check_output(['pkg-config', '--libs', 'x11'], text=True).strip()
    for directory, extra in [
        ('plotlib', ['PLTLIB=libPlt_gDP.a', 'CC=cc', 'CFLAGS=-O2 -DUNDERSCORE ' + cflags]),
        ('bin', ['xfoil', 'FFLOPT=' + flags, 'PLTLIB=' + libraries, 'INSTALLCMD=true']),
    ]:
        with (output / f'{directory}-build.log').open('w') as log:
            subprocess.run(['make', 'FC=' + compiler, 'FFLAGS=' + flags, *extra],
                           cwd=output / 'Xfoil' / directory, stdout=log,
                           stderr=subprocess.STDOUT, check=True)
    (output / 'bin').mkdir(exist_ok=True)
    binary = output / 'bin/xfoil'
    shutil.copy2(output / 'Xfoil/bin/xfoil', binary)
    print(binary)


if __name__ == '__main__':
    build(Path(sys.argv[1]).resolve())
