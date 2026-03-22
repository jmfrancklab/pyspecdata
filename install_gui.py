from PySide6 import QtCore, QtWidgets


def _connect_if_present(signal_owner, name, callback):
    signal = getattr(signal_owner, name, None)
    if signal is not None and hasattr(signal, "connect"):
        signal.connect(callback)


class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(347, 247)
        Dialog.setWindowTitle("Dialog")

        self.gridLayout = QtWidgets.QGridLayout(Dialog)
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setSizeConstraint(QtWidgets.QLayout.SetMaximumSize)

        self.comboBox = QtWidgets.QComboBox(Dialog)
        self.comboBox.setObjectName("comboBox")
        self.comboBox.setEditable(True)
        self.comboBox.addItems(
            [
                "---Select Action to Perform---",
                "Set Data Directory",
                "Set Home Directory",
                "Run Script",
            ]
        )
        self.verticalLayout.addWidget(self.comboBox)

        self.verticalLayout.addItem(
            QtWidgets.QSpacerItem(
                20,
                40,
                QtWidgets.QSizePolicy.Minimum,
                QtWidgets.QSizePolicy.Expanding,
            )
        )

        self.label_3 = QtWidgets.QLabel("<b>Data Directory:</b>", Dialog)
        self.label_3.setObjectName("label_3")
        self.verticalLayout.addWidget(self.label_3)

        self.label_2 = QtWidgets.QLabel("TextLabel", Dialog)
        self.label_2.setObjectName("label_2")
        self.verticalLayout.addWidget(self.label_2)

        self.verticalLayout.addItem(
            QtWidgets.QSpacerItem(
                20,
                40,
                QtWidgets.QSizePolicy.Minimum,
                QtWidgets.QSizePolicy.Expanding,
            )
        )

        self.label_4 = QtWidgets.QLabel("<b>Home Directory:</b>", Dialog)
        self.label_4.setObjectName("label_4")
        self.verticalLayout.addWidget(self.label_4)

        self.label = QtWidgets.QLabel("TextLabel", Dialog)
        self.label.setObjectName("label")
        self.verticalLayout.addWidget(self.label)

        self.verticalLayout.addItem(
            QtWidgets.QSpacerItem(
                20,
                40,
                QtWidgets.QSizePolicy.Minimum,
                QtWidgets.QSizePolicy.Expanding,
            )
        )

        self.gridLayout.addLayout(self.verticalLayout, 0, 0, 1, 1)

        self.buttonBox = QtWidgets.QDialogButtonBox(Dialog)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(
            QtWidgets.QDialogButtonBox.Cancel
            | QtWidgets.QDialogButtonBox.SaveAll
        )
        self.buttonBox.setObjectName("buttonBox")
        self.gridLayout.addWidget(self.buttonBox, 1, 0, 1, 1)

        self.retranslateUi(Dialog)
        self.buttonBox.accepted.connect(Dialog.accept)
        self.buttonBox.rejected.connect(Dialog.reject)
        if hasattr(Dialog, "my_complete"):
            self.buttonBox.accepted.connect(Dialog.my_complete)
        if hasattr(Dialog, "my_change_selection"):
            self.comboBox.currentTextChanged.connect(Dialog.my_change_selection)
        _connect_if_present(Dialog, "my_change_datadir", self.label_2.setText)
        _connect_if_present(Dialog, "my_change_homedir", self.label.setText)
        _connect_if_present(
            Dialog, "my_change_list_selection", self.comboBox.setEditText
        )
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        return
