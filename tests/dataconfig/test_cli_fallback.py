import importlib
import importlib.util
import platform
import sys
from pathlib import Path
from types import ModuleType


ROOT_DIR = Path(__file__).resolve().parents[2]
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

if '_nnls' not in sys.modules:
    sys.modules['_nnls'] = ModuleType('_nnls')


def test_genconfig_without_qt(monkeypatch, tmp_path):
    """Verify that the command falls back to the text file template when Qt is missing."""
    home_dir = tmp_path / "homedir"
    home_dir.mkdir()
    monkeypatch.setenv("HOME", str(home_dir))

    original_import_module = importlib.import_module

    def fake_import_module(name, package=None):
        if name.startswith("PyQt") or name.startswith("PySide"):
            raise ImportError
        return original_import_module(name, package)

    monkeypatch.setattr(importlib, "import_module", fake_import_module)

    package_dir = ROOT_DIR / "pyspecdata"
    if 'pyspecdata' not in sys.modules or not hasattr(sys.modules['pyspecdata'], '__path__'):
        pkg = ModuleType('pyspecdata')
        pkg.__path__ = [str(package_dir)]
        sys.modules['pyspecdata'] = pkg

    def load_module(name, filename):
        spec = importlib.util.spec_from_file_location(name, package_dir / filename)
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
        return module

    load_module('pyspecdata.datadir', 'datadir.py')
    datadir = sys.modules['pyspecdata.datadir']

    datadir.genconfig()

    hide_start = "_" if platform.platform().startswith("Windows") else "."
    config_path = home_dir / f"{hide_start}pyspecdata"
    assert config_path.exists(), "The fallback template should be written when Qt is missing."
    contents = config_path.read_text(encoding="utf-8")
    assert "[General]" in contents
    assert "data_directory =" in contents
    assert "[ExpTypes]" in contents
