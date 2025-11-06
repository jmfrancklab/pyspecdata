import importlib
import platform
from conftest import load_module


def test_genconfig_without_qt(monkeypatch, tmp_path):
    """Verify that the command falls back to the text file template when Qt is
    missing."""
    home_dir = tmp_path / "homedir"
    home_dir.mkdir()
    monkeypatch.setenv("HOME", str(home_dir))

    original_import_module = importlib.import_module

    def fake_import_module(name, package=None):
        if name.startswith("PyQt") or name.startswith("PySide"):
            raise ImportError
        return original_import_module(name, package)

    monkeypatch.setattr(importlib, "import_module", fake_import_module)

    datadir = load_module("datadir")

    datadir.genconfig()

    hide_start = "_" if platform.platform().startswith("Windows") else "."
    config_path = home_dir / f"{hide_start}pyspecdata"
    assert config_path.exists(), (
        "The fallback template should be written when Qt is missing."
    )
    contents = config_path.read_text(encoding="utf-8")
    assert "[General]" in contents
    assert "data_directory =" in contents
    assert "[ExpTypes]" in contents
