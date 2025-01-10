import sys
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QLabel, QLineEdit, QPushButton, QFileDialog, QCheckBox, QHBoxLayout, QMessageBox
import subprocess

class PrimerCheqGUI(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        layout = QVBoxLayout()

        # Primer file
        self.primer_label = QLabel('Primer file:')
        self.primer_input = QLineEdit(self)
        self.primer_button = QPushButton('Browse', self)
        self.primer_button.clicked.connect(self.browse_primer)
        layout.addWidget(self.primer_label)
        layout.addWidget(self.primer_input)
        layout.addWidget(self.primer_button)

        # Prefix
        self.prefix_label = QLabel('Prefix:')
        self.prefix_input = QLineEdit(self)
        layout.addWidget(self.prefix_label)
        layout.addWidget(self.prefix_input)

        # Working directory
        self.working_dir_label = QLabel('Working directory:')
        self.working_dir_input = QLineEdit(self)
        self.working_dir_button = QPushButton('Browse', self)
        self.working_dir_button.clicked.connect(self.browse_working_dir)
        layout.addWidget(self.working_dir_label)
        layout.addWidget(self.working_dir_input)
        layout.addWidget(self.working_dir_button)

        # NCBI Virus
        self.ncbi_virus_label = QLabel('NCBI Virus Taxid:')
        self.ncbi_virus_input = QLineEdit(self)
        layout.addWidget(self.ncbi_virus_label)
        layout.addWidget(self.ncbi_virus_input)

        # NCBI Bacteria
        self.ncbi_bacteria_label = QLabel('NCBI Bacteria Taxid:')
        self.ncbi_bacteria_input = QLineEdit(self)
        layout.addWidget(self.ncbi_bacteria_label)
        layout.addWidget(self.ncbi_bacteria_input)

        # Directory DB
        self.directory_db_label = QLabel('Directory DB:')
        self.directory_db_input = QLineEdit(self)
        self.directory_db_button = QPushButton('Browse', self)
        self.directory_db_button.clicked.connect(self.browse_directory_db)
        layout.addWidget(self.directory_db_label)
        layout.addWidget(self.directory_db_input)
        layout.addWidget(self.directory_db_button)

        # Fasta DB
        self.fasta_db_label = QLabel('Fasta DB:')
        self.fasta_db_input = QLineEdit(self)
        self.fasta_db_button = QPushButton('Browse', self)
        self.fasta_db_button.clicked.connect(self.browse_fasta_db)
        layout.addWidget(self.fasta_db_label)
        layout.addWidget(self.fasta_db_input)
        layout.addWidget(self.fasta_db_button)

        # Glob DB
        self.glob_db_label = QLabel('Glob DB:')
        self.glob_db_input = QLineEdit(self)
        layout.addWidget(self.glob_db_label)
        layout.addWidget(self.glob_db_input)

        # Run button
        self.run_button = QPushButton('Run', self)
        self.run_button.clicked.connect(self.run_script)
        layout.addWidget(self.run_button)

        self.setLayout(layout)
        self.setWindowTitle('Primer Cheq GUI')
        self.show()

    def browse_primer(self):
        fname = QFileDialog.getOpenFileName(self, 'Open file', '', 'FASTA files (*.fasta *.fa)')
        if fname[0]:
            self.primer_input.setText(fname[0])

    def browse_working_dir(self):
        dname = QFileDialog.getExistingDirectory(self, 'Select directory')
        if dname:
            self.working_dir_input.setText(dname)

    def browse_directory_db(self):
        dname = QFileDialog.getExistingDirectory(self, 'Select directory')
        if dname:
            self.directory_db_input.setText(dname)

    def browse_fasta_db(self):
        fname = QFileDialog.getOpenFileName(self, 'Open file', '', 'FASTA files (*.fasta *.fa)')
        if fname[0]:
            self.fasta_db_input.setText(fname[0])

    def run_script(self):
        primer = self.primer_input.text()
        prefix = self.prefix_input.text()
        working_dir = self.working_dir_input.text()
        ncbi_virus = self.ncbi_virus_input.text()
        ncbi_bacteria = self.ncbi_bacteria_input.text()
        directory_db = self.directory_db_input.text()
        fasta_db = self.fasta_db_input.text()
        glob_db = self.glob_db_input.text()

        if not primer or not prefix or not working_dir:
            QMessageBox.warning(self, 'Input Error', 'Please provide all required inputs.')
            return

        cmd = ['python', 'primer_cheq.py', '-p', primer, '-s', prefix, '-w', working_dir]

        if ncbi_virus:
            cmd.extend(['-v', ncbi_virus])
        if ncbi_bacteria:
            cmd.extend(['-b', ncbi_bacteria])
        if directory_db:
            cmd.extend(['-d', directory_db])
        if fasta_db:
            cmd.extend(['-f', fasta_db])
        if glob_db:
            cmd.extend(['-g', glob_db])

        try:
            subprocess.run(cmd, check=True)
            QMessageBox.information(self, 'Success', 'Script ran successfully.')
        except subprocess.CalledProcessError as e:
            QMessageBox.critical(self, 'Error', f'An error occurred: {e}')

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = PrimerCheqGUI()
    sys.exit(app.exec_())