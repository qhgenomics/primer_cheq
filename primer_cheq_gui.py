import sys, os
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QLabel, QLineEdit, QPushButton, QFileDialog, QCheckBox, QGridLayout, QMessageBox
from PyQt5.QtCore import QThread, pyqtSignal, QTimer
import subprocess
from primer_cheq import get_primer_sequences, get_db_glob, get_db_fasta, get_db_fastas, get_db_folder, download_bac, download_virus, blast_primers

class BlastWorker(QThread):
    finished = pyqtSignal()

    def __init__(self, primer_file, primer_dict, working_directory, prefix, blastn_loc):
        super().__init__()
        self.primer_file = primer_file
        self.primer_dict = primer_dict
        self.working_directory = working_directory
        self.prefix = prefix
        self.blastn_loc = blastn_loc

    def run(self):
        blast_primers(self.primer_file, self.primer_dict, self.working_directory, self.prefix, self.blastn_loc)
        self.finished.emit()

class PrimerCheqGUI(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        layout = QGridLayout()
        self.output_label = QLabel('Input options (required)')
        self.output_label.setStyleSheet("background-color: lightyellow; font-weight: bold")
        layout.addWidget(self.output_label, 0, 0, 1, 4)
        # Primer file
        self.primer_label = QLabel('Primer file:')
        self.primer_input = QLineEdit(self)
        self.primer_button = QPushButton('Browse', self)
        self.primer_button.clicked.connect(self.browse_primer)
        layout.addWidget(self.primer_label, 1, 0)
        layout.addWidget(self.primer_input, 1, 1, 1, 2)
        layout.addWidget(self.primer_button, 1, 3)

        self.blast_label = QLabel('Blastn location:')
        self.blast_input = QLineEdit('blastn', self)
        self.blast_button = QPushButton('Browse', self)
        self.blast_button.clicked.connect(self.browse_blast)
        layout.addWidget(self.blast_label, 2, 0)
        layout.addWidget(self.blast_input, 2, 1, 1, 2)
        layout.addWidget(self.blast_button, 2, 3)

        self.output_label = QLabel('Output options (required)')
        self.output_label.setStyleSheet("background-color: lightyellow; font-weight: bold")
        layout.addWidget(self.output_label, 3, 0, 1, 4)
        # Prefix
        self.prefix_label = QLabel('Prefix:')
        self.prefix_input = QLineEdit(self)
        layout.addWidget(self.prefix_label, 4, 0)
        layout.addWidget(self.prefix_input, 4, 1, 1, 3)

        # Working directory
        self.working_dir_label = QLabel('Working directory:')
        self.working_dir_input = QLineEdit(self)
        self.working_dir_button = QPushButton('Browse', self)
        self.working_dir_button.clicked.connect(self.browse_working_dir)
        layout.addWidget(self.working_dir_label, 5, 0)
        layout.addWidget(self.working_dir_input, 5, 1, 1, 2)
        layout.addWidget(self.working_dir_button, 5, 3)

        self.output_label = QLabel('Genome database options (pick 1 or more)')
        self.output_label.setStyleSheet("background-color: lightyellow; font-weight: bold")
        layout.addWidget(self.output_label, 6, 0, 1, 4)
        # NCBI Virus
        self.ncbi_virus_label = QLabel('NCBI Virus Taxid (will download fastas):')
        self.ncbi_virus_input = QLineEdit(self)
        layout.addWidget(self.ncbi_virus_label, 7, 0)
        layout.addWidget(self.ncbi_virus_input, 7, 1, 1, 3)

        # NCBI Bacteria
        self.ncbi_bacteria_label = QLabel('NCBI Bacteria Taxid (will download fastas):')
        self.ncbi_bacteria_input = QLineEdit(self)
        layout.addWidget(self.ncbi_bacteria_label, 8, 0)
        layout.addWidget(self.ncbi_bacteria_input, 8, 1, 1, 3)

        # Directory DB
        self.directory_db_label = QLabel('Directory with multiple fasta files (1 file per genome):')
        self.directory_db_input = QLineEdit(self)
        self.directory_db_button = QPushButton('Browse', self)
        self.directory_db_button.clicked.connect(self.browse_directory_db)
        layout.addWidget(self.directory_db_label, 9, 0)
        layout.addWidget(self.directory_db_input, 9, 1, 1, 2)
        layout.addWidget(self.directory_db_button, 9, 3)

        # Fasta DB
        self.fasta_db_label = QLabel('Fasta file with multiple genomes:')
        self.fasta_db_input = QLineEdit(self)
        self.fasta_db_button = QPushButton('Browse', self)
        self.fasta_db_button.clicked.connect(self.browse_fasta_db)
        layout.addWidget(self.fasta_db_label, 10, 0)
        layout.addWidget(self.fasta_db_input, 10, 1, 1, 2)
        layout.addWidget(self.fasta_db_button, 10, 3)

        # Glob DB
        self.glob_db_label = QLabel('Path with wildcard (*) that resolves to one or more genomes:')
        self.glob_db_input = QLineEdit(self)
        layout.addWidget(self.glob_db_label, 11, 0)
        layout.addWidget(self.glob_db_input, 11, 1, 1, 3)

        # Run button
        self.run_button = QPushButton('Run', self)
        self.run_button.clicked.connect(self.run_script)
        layout.addWidget(self.run_button, 12, 0, 1, 4)

        self.setLayout(layout)
        self.setWindowTitle('Primer Cheq GUI')
        self.show()

    def browse_primer(self):
        fname = QFileDialog.getOpenFileName(self, 'Open file', '', 'FASTA files (*.fasta *.fa)')
        if fname[0]:
            self.primer_input.setText(fname[0])

    def browse_blast(self):
        fname = QFileDialog.getOpenFileName(self, 'Open file', '', 'BLASTn binary (*.exe, *)')
        if fname[0]:
            self.blast_input.setText(fname[0])

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
        primers = self.primer_input.text().strip()
        prefix = self.prefix_input.text().strip()
        working_directory = self.working_dir_input.text().strip()
        blastn_loc = '"' + self.blast_input.text().strip() + '"'
        ncbi_virus = self.ncbi_virus_input.text().strip()
        ncbi_bacteria = self.ncbi_bacteria_input.text().strip()
        directory_db = self.directory_db_input.text().strip()
        fasta_db = self.fasta_db_input.text().strip()
        glob_db = self.glob_db_input.text().strip()

        if not primers or not prefix or not working_directory:
            QMessageBox.warning(self, 'Input Error', 'Please provide all required inputs.')
            return
        if not ncbi_virus and not ncbi_bacteria and not directory_db and not fasta_db and not glob_db:
            QMessageBox.warning(self, 'Input Error', 'Please provide at least one genome database.')

        primer_dict = get_primer_sequences(primers)
        if not os.path.exists(working_directory):
            os.makedirs(working_directory)
        elif os.path.exists(working_directory) and not os.path.isdir(working_directory):
            QMessageBox.critical(self, 'Error', 'An error occurred: working directory is not a directory.')

        primer_file = os.path.join(working_directory, prefix + "_db.fasta")
        open(primer_file, 'w').close()
        if ncbi_virus:
            fasta_file = download_virus(ncbi_virus, working_directory, prefix)
            get_db_fasta(fasta_file, working_directory, prefix)
        if ncbi_bacteria:
            fasta_files = download_bac(ncbi_bacteria, working_directory, prefix)
            get_db_fastas(fasta_files, working_directory, prefix)
        if directory_db:
            fasta_files = get_db_folder(directory_db)
            get_db_fastas(fasta_files, working_directory, prefix)
        if fasta_db:
            get_db_fasta(fasta_db, working_directory, prefix)
        if glob_db:
            fasta_files = get_db_glob(glob_db)
            get_db_fastas(fasta_files, working_directory, prefix)

        self.worker = BlastWorker(primer_file, primer_dict, working_directory, prefix, blastn_loc)
        self.worker.finished.connect(self.on_blast_finished)
        self.worker.start()

        self.animation_label = QLabel('Running...')
        self.animation_label.setStyleSheet("background-color: lightgreen; font-weight: bold")
        self.layout().addWidget(self.animation_label, 13, 0, 1, 4)
        self.timer = QTimer(self)
        self.timer.timeout.connect(self.update_animation)
        self.timer.start(500)

    def update_animation(self):
        current_text = self.animation_label.text()
        if current_text.endswith('...'):
            self.animation_label.setText('Running')
        else:
            self.animation_label.setText(current_text + '.')

    def on_blast_finished(self):
        self.timer.stop()
        self.animation_label.setText('Finished')
        QMessageBox.information(self, 'Success', 'Script ran successfully.')

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = PrimerCheqGUI()
    sys.exit(app.exec_())