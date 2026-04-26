# -*- mode: python ; coding: utf-8 -*-
"""MultiGaussFit — onefile spec (portable, slower startup)"""

EXCLUDES = [
    'PyQt5', 'PyQt6', 'PySide2', 'PySide6', 'wx', 'gtk', 'gi',
    'matplotlib.backends.backend_qt5agg', 'matplotlib.backends.backend_qt5',
    'matplotlib.backends.backend_qt', 'matplotlib.backends.backend_qtagg',
    'matplotlib.backends.backend_webagg', 'matplotlib.backends.backend_nbagg',
    'matplotlib.backends.backend_gtk3agg', 'matplotlib.backends.backend_gtk3',
    'matplotlib.backends.backend_gtk4agg', 'matplotlib.backends.backend_gtk4',
    'matplotlib.backends.backend_wx', 'matplotlib.backends.backend_wxagg',
    'matplotlib.backends.backend_cairo', 'matplotlib.backends.backend_macosx',
    'IPython', 'jupyter', 'notebook', 'nbconvert', 'nbformat',
    'docutils', 'sphinx', 'cv2', 'skimage', 'sklearn',
    'sqlalchemy', 'tornado',
]

a = Analysis(['MultiGaussFit.py'], excludes=EXCLUDES,
             hiddenimports=['matplotlib.backends.backend_tkagg'],
             datas=[('logo.ico', '.'), ('logo.png', '.')])
a.binaries = [b for b in a.binaries if not b[0].endswith(('.pdb', '.lib'))]
pyz = PYZ(a.pure)
exe = EXE(pyz, a.scripts, a.binaries, a.datas, [],
          name='MultiGaussFit', strip=False, upx=False, console=False, icon='logo.ico', version='file_version_info.txt')
