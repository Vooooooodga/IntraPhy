#!/usr/bin/env python3
"""Check the installed distribution from outside the source tree."""
import importlib.metadata
import importlib.util
import shutil
import subprocess
import sys


def main():
    import intraphy
    distribution=importlib.metadata.distribution('intraphy')
    entries=[entry.name for entry in distribution.entry_points if entry.group=='console_scripts']
    assert entries==['intraphy'], entries
    assert intraphy.__version__ == distribution.version, (intraphy.__version__, distribution.version)
    previous='insi'+'phy'
    assert importlib.util.find_spec(previous) is None, 'A separately installed old package remains in the environment'
    assert shutil.which(previous) is None, 'A separately installed old executable remains on PATH'
    for command in ([sys.executable,'-I','-m','intraphy','--version'],['intraphy','--version']):
        completed=subprocess.run(command,check=True,capture_output=True,text=True)
        assert completed.stdout.strip() == f'intraphy {distribution.version}', completed.stdout
    print('Installed namespace and entry points: passed')

if __name__=='__main__':
    main()
