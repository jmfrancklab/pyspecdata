import configparser
import importlib
import importlib.util
import platform
import subprocess
import sys
from pathlib import Path
from types import ModuleType, SimpleNamespace

import pytest


ROOT_DIR = Path(__file__).resolve().parents[2]
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

if '_nnls' not in sys.modules:
    sys.modules['_nnls'] = ModuleType('_nnls')


def test_genconfig_with_qt(monkeypatch, tmp_path):
    """Exercise the Qt configuration editor and verify the resulting file."""
    home_dir = tmp_path / "homedir"
    home_dir.mkdir()
    sample_path = Path(__file__).parent / "sample.pyspecdata"
    hide_start = "_" if platform.platform().startswith("Windows") else "."
    config_path = home_dir / f"{hide_start}pyspecdata"
    config_path.write_text(sample_path.read_text(encoding="utf-8"), encoding="utf-8")
    monkeypatch.setenv("HOME", str(home_dir))
    monkeypatch.setenv("QT_QPA_PLATFORM", "offscreen")

    original_check_output = subprocess.check_output

    def fake_check_output(cmd, *args, **kwargs):
        if isinstance(cmd, list) and len(cmd) >= 2 and cmd[0] == "rclone" and cmd[1] == "listremotes":
            return "jmf_teams:\ncornell_box:\n".encode("utf-8")
        return original_check_output(cmd, *args, **kwargs)

    monkeypatch.setattr(subprocess, "check_output", fake_check_output)

    class StubSignal:
        def __init__(self):
            self._callbacks = []

        def connect(self, callback):
            self._callbacks.append(callback)

        def emit(self, *args, **kwargs):
            for callback in list(self._callbacks):
                callback(*args, **kwargs)

    class StubStyle:
        def standardIcon(self, *args, **kwargs):
            return None

    class StubStyleClass:
        SP_DirOpenIcon = 0

    class BaseWidget:
        def __init__(self, *args, **kwargs):
            self._style = StubStyle()

        def style(self):
            return self._style

        def setLayout(self, layout):
            self.layout = layout

        def setSizePolicy(self, *args, **kwargs):
            return

    class StubWidget(BaseWidget):
        pass

    class StubDialog(StubWidget):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self._accepted = False

        def setWindowTitle(self, title):
            self.window_title = title

        def resize(self, *args, **kwargs):
            self.size = (args, kwargs)

        def accept(self):
            self._accepted = True

        def reject(self):
            self._accepted = False

        def exec_(self):
            return 1 if self._accepted else 0

    class StubApplication:
        _instance = None

        def __init__(self, *args, **kwargs):
            StubApplication._instance = self

        @classmethod
        def instance(cls):
            return cls._instance

        def exec_(self):
            return 0

    class StubLayout:
        def __init__(self, parent=None):
            self.items = []
            if parent is not None:
                parent.layout = self

        def addWidget(self, widget):
            self.items.append(widget)

        def addLayout(self, layout):
            self.items.append(layout)

        def addStretch(self, *args, **kwargs):
            return

        def insertWidget(self, index, widget):
            if index < 0 or index > len(self.items):
                self.items.append(widget)
            else:
                self.items.insert(index, widget)

        def setContentsMargins(self, *args, **kwargs):
            return

        def setSpacing(self, *args, **kwargs):
            return

        def count(self):
            return len(self.items)

    class StubVBoxLayout(StubLayout):
        pass

    class StubHBoxLayout(StubLayout):
        pass

    class StubTabWidget(StubWidget):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.tabs = []

        def addTab(self, widget, name):
            self.tabs.append((widget, name))

    class StubScrollArea(StubWidget):
        def setWidgetResizable(self, *args, **kwargs):
            return

        def setWidget(self, widget):
            self.widget = widget

    class StubLineEdit(StubWidget):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self._text = ''
            self.enabled = True

        def setText(self, value):
            self._text = str(value)

        def text(self):
            return self._text

        def clear(self):
            self._text = ''

        def setEnabled(self, state):
            self.enabled = state

    class StubLabel(StubWidget):
        def __init__(self, text='', *args, **kwargs):
            super().__init__(*args, **kwargs)
            self._text = text

    class StubButton(StubWidget):
        def __init__(self, text='', *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.text = text
            self.clicked = StubSignal()

        def click(self):
            self.clicked.emit()

        def setIcon(self, *args, **kwargs):
            return

    class StubComboBox(StubWidget):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.items = []
            self.index = 0
            self.currentIndexChanged = StubSignal()

        def addItem(self, text):
            self.items.append(text)

        def setCurrentIndex(self, index):
            if 0 <= index < len(self.items):
                self.index = index
            else:
                self.index = 0 if self.items else 0
            self.currentIndexChanged.emit(self.index)

        def currentText(self):
            if not self.items:
                return ''
            return self.items[self.index]

        def findText(self, text):
            try:
                return self.items.index(text)
            except ValueError:
                return -1

    class StubFileDialog:
        @staticmethod
        def getExistingDirectory(*args, **kwargs):
            return ''

    class StubSizePolicy:
        Minimum = 0
        Expanding = 1

    stub_widgets = SimpleNamespace(
        QApplication=StubApplication,
        QWidget=StubWidget,
        QDialog=StubDialog,
        QVBoxLayout=StubVBoxLayout,
        QHBoxLayout=StubHBoxLayout,
        QTabWidget=StubTabWidget,
        QScrollArea=StubScrollArea,
        QLineEdit=StubLineEdit,
        QLabel=StubLabel,
        QPushButton=StubButton,
        QToolButton=StubButton,
        QComboBox=StubComboBox,
        QFileDialog=StubFileDialog,
        QSizePolicy=StubSizePolicy,
        QStyle=StubStyleClass,
    )

    original_import_module = importlib.import_module

    def fake_import_module(name, package=None):
        if name == 'PyQt5.QtWidgets':
            return stub_widgets
        return original_import_module(name, package)

    monkeypatch.setattr(importlib, 'import_module', fake_import_module)

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
    latexscripts = load_module('pyspecdata.latexscripts', 'latexscripts.py')

    def fake_exec(self):
        for row in self.file_rows:
            if row['name_edit'].text() == 'ag_processed_data':
                index = row['remote_combo'].findText('cornell_box')
                assert index != -1
                row['remote_combo'].setCurrentIndex(index)
                row['remote_path_edit'].setText('exp_data/AG_processed_data')
            if row['name_edit'].text() == 'francklab_esr/romana':
                index = row['remote_combo'].findText('None')
                assert index != -1
                row['remote_combo'].setCurrentIndex(index)
        self.add_button.click()
        new_row = self.file_rows[-1]
        new_row['name_edit'].setText('new_entry')
        new_row['path_edit'].setText('/tmp/new_entry')
        for entry in self.general_entries:
            if entry['key_edit'].text() == 'qesr conversion':
                entry['value_edit'].setText('999')
        self.data_directory_edit.setText('/tmp/data_dir')
        self.save_button.click()
        return 1

    monkeypatch.setattr(stub_widgets.QDialog, 'exec_', fake_exec)

    latexscripts.genconfig()

    config = configparser.ConfigParser()
    config.read(config_path, encoding='utf-8')

    assert config.get('General', 'data_directory') == '/tmp/data_dir'
    assert config.get('General', 'qesr conversion') == '999'
    assert config.get('ExpTypes', 'new_entry') == '/tmp/new_entry'
    assert config.get('RcloneRemotes', 'ag_processed_data') == 'cornell_box:exp_data/AG_processed_data'
    assert not config.has_option('RcloneRemotes', 'francklab_esr/romana')
    assert config.get('RcloneRemotes', '2deldor') == 'cornell_box:exp_data/2DELDOR'
    assert config.has_section('zenodo'), 'Existing sections outside the editor must remain.'
